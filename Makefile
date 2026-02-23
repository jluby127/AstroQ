# Simple AstroQ Makefile
# Creates directory structure and copies config template with current date

# Configuration variables
# Conda environment to use for all commands
CONDA_ENV ?= astroq
WORKDIR ?= $(CC_OUTPUT_PATH)
HOLDERS_DIR ?= $(CC_RESULTS_COPY_PATH)/$(SEMESTER)/$(DATE)
SEMESTER_DIR = $(WORKDIR)/$(SEMESTER)
DATE_DIR = $(SEMESTER_DIR)/$(DATE)

# Changed only once per semester
SEMESTER ?= 2026A
START_DATE ?= 2026-02-01
END_DATE ?= 2026-07-31
BANDS ?= band1 band2 band3 full-band1 full-band2 full-band3
FILLER_PROGRAM ?= 2026A_E475

# Cross-platform sed in-place flag (Darwin uses -i '')
UNAME_S := $(shell uname)
SED_I_FLAG :=
ifeq ($(UNAME_S),Darwin)
  SED_I_FLAG := ''
endif

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
	@echo "📋 Copying Magiq script files to holders..."
	@$(MAKE) copy_magiq_scripts
	@echo "📝 Combining all band logs..."
	@$(MAKE) combine_logs
	@echo "🔍 Checking night plans..."
	@$(MAKE) check_night_plans

# Final target for each band - depends on plan-night completion
$(DATE_DIR)/%/plan-night-complete: $(DATE_DIR)/%/plan-night-run
	@echo "✅ Band $(notdir $(@D)) workflow complete!"
	@touch $@

# Run plan-night command
$(DATE_DIR)/%/plan-night-run: $(DATE_DIR)/%/plan-semester-run
	@echo "🌙 Running plan-night for band $(notdir $(@D))..."
	@cd $(@D) && conda run -n $(CONDA_ENV) astroq plan-night -cf config.ini 2>&1 | tee -a astroq.log
	@touch $@

# Run plan-semester command (unified for all bands)
$(DATE_DIR)/%/plan-semester-run: $(DATE_DIR)/%/prep-run
	@echo "📅 Running plan-semester for band $(notdir $(@D))..."
	@cd $(@D) && conda run -n $(CONDA_ENV) astroq plan-semester -cf config.ini 2>&1 | tee -a astroq.log
	@touch $@

# Unified prep command for all bands
$(DATE_DIR)/%/prep-run: $(DATE_DIR)/%/config.ini
	@echo "🔧 Running prep for band $(notdir $(@D))..."
	@BAND_NUM=$$(echo $(notdir $(@D)) | sed 's/band//' | sed 's/full-//'); \
	IS_FULL_BAND=$$(echo $(notdir $(@D)) | grep -q '^full-' && echo "true" || echo "false"); \
	if [ "$$IS_FULL_BAND" = "true" ]; then \
		echo "📅 Running prep for full-band $$BAND_NUM..." && \
		cd $(@D) && conda run -n $(CONDA_ENV) astroq prep kpfcc -cf config.ini -fillers $(FILLER_PROGRAM) -band $$BAND_NUM -full 2>&1 | tee -a astroq.log; \
	else \
		echo "📊 Running prep for band $$BAND_NUM..." && \
		cd $(@D) && conda run -n $(CONDA_ENV) astroq prep kpfcc -cf config.ini -fillers $(FILLER_PROGRAM) -band $$BAND_NUM 2>&1 | tee -a astroq.log; \
	fi
	@touch $@

# Create config file for a specific band
$(DATE_DIR)/%/config.ini: create_dirs
	@echo "📋 Copying config template for band $(notdir $(@D))..."
	@mkdir -p $(@D)
	@cp config_template.ini $@
	@sed -i $(SED_I_FLAG) "s|CURRENT_DATE_PLACEHOLDER|$(DATE)|g" $@
	@sed -i $(SED_I_FLAG) "s|START_DATE_PLACEHOLDER|$(START_DATE)|g" $@
	@sed -i $(SED_I_FLAG) "s|END_DATE_PLACEHOLDER|$(END_DATE)|g" $@
	@sed -i $(SED_I_FLAG) "s|SEMESTER_PLACEHOLDER|$(SEMESTER)|g" $@
	@sed -i $(SED_I_FLAG) "s|WORKDIR_PLACEHOLDER|$(@D)|g" $@
	# If this is band3 or full-band3, tighten the max_solve_gap for semester solver
	@BAND_NAME=$(notdir $(@D)); \
	if [ "$$BAND_NAME" = "band3" ] || [ "$$BAND_NAME" = "full-band3" ]; then \
		echo "⚙️  Detected $$BAND_NAME: setting [semester] max_solve_gap = 0.0005"; \
		sed -i $(SED_I_FLAG) -e '/^\[semester\]/,/^\[/{s/^max_solve_gap.*/max_solve_gap = 0.0005/;}' $@; \
	fi
	# If this is a full band (full-band1, full-band2, full-band3), set max_solve_time for night solver
	@BAND_NAME=$(notdir $(@D)); \
	if echo "$$BAND_NAME" | grep -q "^full-band"; then \
		echo "⚙️  Detected $$BAND_NAME: setting [night] max_solve_time = 300"; \
		sed -i $(SED_I_FLAG) -e '/^\[night\]/,/^\[/{s/^max_solve_time.*/max_solve_time = 300/;}' $@; \
	fi
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

# Copy script_{date}_nominal.txt to holders as magiq_backup_{band}.txt
copy_magiq_scripts:
	@echo "📋 Copying Magiq script files to holders directories..."
	@for band in $(BANDS); do \
		echo "📋 Copying script for band $$band..."; \
		mkdir -p $(HOLDERS_DIR)/$$band/output; \
		cp $(DATE_DIR)/$$band/outputs/script_$(DATE)_nominal.txt $(HOLDERS_DIR)/$$band/output/magiq_backup_$$band.txt; \
		echo "✅ Magiq script copied for band $$band"; \
	done
	@echo "✅ All Magiq script files copied to holders directories!"

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

# Check night plans
check_night_plans:
	@echo "🔍 Running check_night_plans.py..."
	@cd $(shell dirname $(firstword $(MAKEFILE_LIST))) && conda run -n $(CONDA_ENV) python check_night_plans.py -s $(SEMESTER) -d $(DATE) $(HOLDERS_DIR)
	@echo "✅ Night plans check complete!"

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

.PHONY: all create_dirs clean status copy_observe_orders copy_only webapp complete combine_logs check_night_plans 