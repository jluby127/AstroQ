# Simple AstroQ Makefile
# Creates directory structure and copies config template with current date

# Configuration variables
DATE ?= 2025-08-28
# Changed only once per semester
SEMESTER ?= 2025B
START_DATE ?= 2025-08-01
END_DATE ?= 2026-01-31
BANDS ?= band1 band3
WORKDIR ?= /Users/jack/Desktop
# WORKDIR ?= /home/kpfcc/AstroQ
RUN_SCRIPT_PATH ?= /path/to/run.sh

# Derived paths
SEMESTER_DIR = $(WORKDIR)/$(SEMESTER)
DATE_DIR = $(SEMESTER_DIR)/$(DATE)
HOLDERS_DIR = $(WORKDIR)/holders/$(SEMESTER)/$(DATE)
# HOLDERS_DIR = /s/sdata1701/Schedules/$(SEMESTER)/$(DATE)

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
	@$(MAKE) copy_custom_csv

# Final target for each band - depends on plan-night completion
$(DATE_DIR)/%/plan-night-complete: $(DATE_DIR)/%/plan-night-run
	@echo "✅ Band $(notdir $(@D)) workflow complete!"
	@touch $@

# Run plan-night command
$(DATE_DIR)/%/plan-night-run: $(DATE_DIR)/%/plan-semester-run
	@echo "🌙 Running plan-night for band $(notdir $(@D))..."
	@cd $(@D) && conda run -n astroq astroq kpfcc plan-night -cf config.ini
	@touch $@

# Run plan-semester command for band3 (specific target)
$(DATE_DIR)/band3/plan-semester-run: $(DATE_DIR)/band3/prep-run
	@echo "📅 Running plan-semester for band3 with -b3 True..."
	@cd $(@D) && conda run -n astroq astroq kpfcc plan-semester -cf config.ini -b3 True
	@touch $@

# Run plan-semester command for band1 (specific target)
$(DATE_DIR)/band1/plan-semester-run: $(DATE_DIR)/band1/prep-run
	@echo "📅 Running plan-semester for band1 without -b3 flag..."
	@cd $(@D) && conda run -n astroq astroq kpfcc plan-semester -cf config.ini
	@touch $@

# Run prep command
$(DATE_DIR)/%/prep-run: $(DATE_DIR)/%/config.ini
	@echo "🔧 Running prep for band $(notdir $(@D))..."
	@cd $(@D) && conda run -n astroq astroq kpfcc prep -cf config.ini
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

# Create directory structure
create_dirs:
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

# Copy WORKDIR/custom.csv to {semester}/{date}/{band}/custom.csv for all bands except band3
copy_custom_csv:
	@for band in $(BANDS); do \
		if [ "$$band" != "band3" ]; then \
			echo "Copying custom.csv to $(DATE_DIR)/$$band/"; \
			cp $(SEMESTER_DIR)/custom.csv $(DATE_DIR)/$$band/custom.csv; \
		fi \
	done

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
	@echo "  Date: $(DATE)"
	@echo "  Start Date: $(START_DATE)"
	@echo "  End Date: $(END_DATE)"
	@echo "  Bands: $(BANDS)"
	@echo "  Work Directory: $(WORKDIR)"
	@echo "  Date Directory: $(DATE_DIR)"

# Complete workflow (workflow + copy + webapp)
complete: all
	@echo "📋 Creating holders directories and copying ObserveOrder files..."
	@$(MAKE) copy_observe_orders
	@echo "🌐 Launching webapp..."
	@$(MAKE) webapp

.PHONY: all create_dirs clean status copy_observe_orders copy_only webapp complete copy_custom_csv 
