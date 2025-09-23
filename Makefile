# Simple AstroQ Makefile
# Creates directory structure and copies config template with current date

# Configuration variables
# Conda environment to use for all commands
CONDA_ENV ?= astroq_new
WORKDIR ?= $(CC_OUTPUT_PATH)
HOLDERS_DIR ?= $(CC_OUTPUT_PATH)/$(CC_RESULTS_COPY_PATH)/$(SEMESTER)/$(DATE)
SEMESTER_DIR = $(WORKDIR)/$(SEMESTER)
DATE_DIR = $(SEMESTER_DIR)/$(DATE)

# Changed only once per semester
SEMESTER ?= 2025B
START_DATE ?= 2025-08-01
END_DATE ?= 2026-01-31
BANDS ?= full-band1 #band1 band2 band3 full-band1 full-band2 full-band3
FILLER_PROGRAM ?= 2025B_E473

# Date validation and override
DATE ?= $(shell date +%Y-%m-%d)
ifdef date
	DATE := $(date)
endif

# Mark config.ini files as precious to prevent automatic deletion
.PRECIOUS: $(foreach band,$(BANDS),$(DATE_DIR)/$(band)/config.ini)

# Default target - run all bands through the complete workflow
all: $(foreach band,$(BANDS),$(DATE_DIR)/$(band)/plan-night-complete)
	@echo "✅ Complete workflow finished for all bands!"
	@echo "📁 Results available in:"
	@for band in $(BANDS); do \
		echo "  $(DATE_DIR)/$$band/"; \
	done
	@echo "📋 Creating holders directories and copying ObserveOrder files..."
	@$(MAKE) copy_observe_orders
	@echo "📝 Combining all band logs..."
	@$(MAKE) combine_logs

# Final target for each band - depends on plan-night completion
$(DATE_DIR)/%/plan-night-complete: $(DATE_DIR)/%/plan-night-run
	@echo "✅ Band $(notdir $(@D)) workflow complete!"
	@touch $@

# Run plan-night command
$(DATE_DIR)/%/plan-night-run: $(DATE_DIR)/%/plan-semester-run
	@echo "🌙 Running plan-night for band $(notdir $(@D))..."
	@cd $(@D) && conda run -n $(CONDA_ENV) astroq kpfcc plan-night -cf config.ini 2>&1 | tee -a astroq.log
	@touch $@

# Run plan-semester command (unified for all bands)
$(DATE_DIR)/%/plan-semester-run: $(DATE_DIR)/%/prep-run
	@echo "📅 Running plan-semester for band $(notdir $(@D))..."
	@cd $(@D) && conda run -n $(CONDA_ENV) astroq kpfcc plan-semester -cf config.ini 2>&1 | tee -a astroq.log
	@touch $@

# Unified prep command for all bands
$(DATE_DIR)/%/prep-run: $(DATE_DIR)/%/config.ini
	@echo "🔧 Running prep for band $(notdir $(@D))..."
	@BAND_NUM=$$(echo $(notdir $(@D)) | sed 's/band//' | sed 's/full-//'); \
	IS_FULL_BAND=$$(echo $(notdir $(@D)) | grep -q '^full-' && echo "true" || echo "false"); \
	if [ "$$IS_FULL_BAND" = "true" ]; then \
		echo "📅 Running prep for full-band $$BAND_NUM..." && \
		cd $(@D) && conda run -n $(CONDA_ENV) astroq kpfcc prep -cf config.ini -fillers $(FILLER_PROGRAM) -band $$BAND_NUM -full 2>&1 | tee -a astroq.log; \
	else \
		echo "📊 Running prep for band $$BAND_NUM..." && \
		cd $(@D) && conda run -n $(CONDA_ENV) astroq kpfcc prep -cf config.ini -fillers $(FILLER_PROGRAM) -band $$BAND_NUM 2>&1 | tee -a astroq.log; \
	fi
	@touch $@

# Create config file for a specific band
$(DATE_DIR)/%/config.ini: create_dirs
	@echo "📋 Copying config template for band $(notdir $(@D))..."
	@mkdir -p $(@D)
	@cp config_template.ini $@
	@echo "📝 Updating placeholders for band $(notdir $(@D))..."
	# @sed -i "s|CURRENT_DATE_PLACEHOLDER|$(DATE)|g" $@
	# @sed -i "s|START_DATE_PLACEHOLDER|$(START_DATE)|g" $@
	# @sed -i "s|END_DATE_PLACEHOLDER|$(END_DATE)|g" $@
	# @sed -i "s|SEMESTER_PLACEHOLDER|$(SEMESTER)|g" $@
	# @sed -i "s|WORKDIR_PLACEHOLDER|$(@D)|g" $@
	@sed -i '' "s|CURRENT_DATE_PLACEHOLDER|$(DATE)|g" $@
	@sed -i '' "s|START_DATE_PLACEHOLDER|$(START_DATE)|g" $@
	@sed -i '' "s|END_DATE_PLACEHOLDER|$(END_DATE)|g" $@
	@sed -i '' "s|SEMESTER_PLACEHOLDER|$(SEMESTER)|g" $@
	@sed -i '' "s|WORKDIR_PLACEHOLDER|$(@D)|g" $@
	@echo "✅ Config file created and updated for band $(notdir $(@D))"

# Validate date format
validate_date:
ifdef date
	@echo "🔍 Validating custom date: $(date)"
	@if ! echo "$(date)" | grep -E '^[0-9]{4}-[0-9]{2}-[0-9]{2}$$' > /dev/null; then \
		echo "❌ Error: Invalid date format: $(date)"; \
		echo "   Please use YYYY-MM-DD format (e.g., 2025-09-16)"; \
		exit 1; \
	fi
	@if ! python3 -c "from datetime import datetime; datetime.strptime('$(date)', '%Y-%m-%d')" > /dev/null 2>&1; then \
		echo "❌ Error: Invalid date: $(date)"; \
		echo "   Please provide a valid date in YYYY-MM-DD format"; \
		exit 1; \
	fi
	@echo "✅ Date validation passed: $(date)"
else
	@echo "📅 Using current date: $(DATE)"
endif

# Create directory structure
create_dirs: validate_date
	@echo "📁 Creating directory structure..."
	@mkdir -p $(DATE_DIR)
	@echo "✅ Directories created"

# Copy ObserveOrder files to holders directories (doesn't depend on workflow)
copy_observe_orders:
	@echo "📋 Copying ObserveOrder files to holders directories..."
	@for band in $(BANDS); do \
		echo "📋 Copying ObserveOrder file for band $$band..."; \
		mkdir -p $(HOLDERS_DIR)/$$band/output; \
		cp $(DATE_DIR)/$$band/outputs/ObserveOrder_$(DATE).txt $(HOLDERS_DIR)/$$band/output/night_plan.csv; \
		echo "✅ ObserveOrder file copied for band $$band"; \
	done
	@echo "✅ All ObserveOrder files copied to holders directories!"

# Combine all band logs into a single file
combine_logs:
	@echo "📝 Combining logs from all bands..."
	@echo "# Combined AstroQ Logs for $(DATE)" > $(DATE_DIR)/all_band_logs.log
	@echo "# Generated on: $$(date)" >> $(DATE_DIR)/all_band_logs.log
	@echo "# Bands: $(BANDS)" >> $(DATE_DIR)/all_band_logs.log
	@echo "" >> $(DATE_DIR)/all_band_logs.log
	@for band in $(BANDS); do \
		if [ -f "$(DATE_DIR)/$$band/astroq.log" ]; then \
			echo "## Band: $$band" >> $(DATE_DIR)/all_band_logs.log; \
			echo "## Log file: $(DATE_DIR)/$$band/astroq.log" >> $(DATE_DIR)/all_band_logs.log; \
			echo "## Timestamp: $$(date)" >> $(DATE_DIR)/all_band_logs.log; \
			echo "" >> $(DATE_DIR)/all_band_logs.log; \
			cat $(DATE_DIR)/$$band/astroq.log >> $(DATE_DIR)/all_band_logs.log; \
			echo "" >> $(DATE_DIR)/all_band_logs.log; \
			echo "## End of $$band log" >> $(DATE_DIR)/all_band_logs.log; \
			echo "" >> $(DATE_DIR)/all_band_logs.log; \
			echo "----------------------------------------" >> $(DATE_DIR)/all_band_logs.log; \
			echo "" >> $(DATE_DIR)/all_band_logs.log; \
		else \
			echo "⚠️  Warning: No log file found for band $$band at $(DATE_DIR)/$$band/astroq.log"; \
		fi \
	done
	@echo "✅ Combined log file created: $(DATE_DIR)/all_band_logs.log"

# Launch webapp
webapp:
	@echo "🌐 Launching AstroQ webapp..."
	@echo "📁 Using workdir as uptree path: $(WORKDIR)"
	@echo "🚀 Running launch_webapp.sh from $(WORKDIR)..."
	@chmod +x $(WORKDIR)/launch_webapp.sh
	@$(WORKDIR)/launch_webapp.sh

# Clean up
clean:
	@echo "🧹 Cleaning up..."
	@echo "⚠️  This will remove: $(DATE_DIR)"
	@echo -n "Are you sure? (y/N): " && read -r confirm && [ "$$confirm" = "y" ] || exit 1
	@rm -rf $(DATE_DIR)
	@echo "✅ Cleanup complete"

# Show status
status:
	@echo "📊 Status:"
	@echo "  Semester: $(SEMESTER)"
	@if [ -n "$(date)" ]; then \
		echo "  Date: $(DATE) (custom)"; \
	else \
		echo "  Date: $(DATE) (current)"; \
	fi
	@echo "  Start Date: $(START_DATE)"
	@echo "  End Date: $(END_DATE)"
	@echo "  Bands: $(BANDS)"
	@echo "  Filler Program: $(FILLER_PROGRAM)"
	@echo "  Work Directory: $(WORKDIR)"
	@echo "  Date Directory: $(DATE_DIR)"

# Complete workflow (workflow + copy + webapp)
complete: all
	@echo "📋 Creating holders directories and copying ObserveOrder files..."
	@$(MAKE) copy_observe_orders
	@echo "🌐 Launching webapp..."
	@$(MAKE) webapp

.PHONY: all create_dirs clean status copy_observe_orders copy_only webapp complete combine_logs 