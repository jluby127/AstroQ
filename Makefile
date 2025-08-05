# Simple AstroQ Makefile
# Creates directory structure and copies config template with current date

# Configuration variables
DATE ?= 2025-08-02
# Changed only once per semester
SEMESTER ?= 2025B
START_DATE ?= 2025-08-01
END_DATE ?= 2026-01-31
BANDS ?= band1 band3
WORKDIR ?= /Users/jack/Desktop

# Derived paths
SEMESTER_DIR = $(WORKDIR)/$(SEMESTER)
DATE_DIR = $(SEMESTER_DIR)/$(DATE)
HOLDERS_DIR = $(WORKDIR)/holders/$(SEMESTER)/$(DATE)

# Default target - run all bands through the complete workflow
all: $(foreach band,$(BANDS),$(DATE_DIR)/$(band)/plan-night-complete)
	@echo "✅ Complete workflow finished for all bands!"
	@echo "📁 Results available in:"
	@for band in $(BANDS); do \
		echo "  $(DATE_DIR)/$$band/"; \
	done
	@echo "📋 Creating holders directories and copying ObserveOrder files..."
	@$(MAKE) copy_observe_orders

# Final target for each band - depends on plan-night completion
$(DATE_DIR)/%/plan-night-complete: $(DATE_DIR)/%/plan-night-run
	@echo "✅ Band $(notdir $(@D)) workflow complete!"
	@touch $@

# Run plan-night command
$(DATE_DIR)/%/plan-night-run: $(DATE_DIR)/%/plan-semester-run
	@echo "🌙 Running plan-night for band $(notdir $(@D))..."
	@cd $(@D) && conda run -n astroq astroq kpfcc plan-night -cf config.ini
	@touch $@

# Run plan-semester command
$(DATE_DIR)/%/plan-semester-run: $(DATE_DIR)/%/prep-run
	@echo "📅 Running plan-semester for band $(notdir $(@D))..."
	@cd $(@D) && if [ "$(notdir $(@D))" = "band3" ]; then \
		conda run -n astroq astroq kpfcc plan-semester -cf config.ini -b3 True; \
	else \
		conda run -n astroq astroq kpfcc plan-semester -cf config.ini; \
	fi
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

# Copy ObserveOrder files to holders directories
copy_observe_orders: $(foreach band,$(BANDS),$(HOLDERS_DIR)/$(band)/outputs/night_plan.csv)
	@echo "✅ All ObserveOrder files copied to holders directories!"

# Copy ObserveOrder file for a specific band
$(HOLDERS_DIR)/%/outputs/night_plan.csv: $(DATE_DIR)/%/plan-night-complete
	@echo "📋 Copying ObserveOrder file for band $(notdir $(patsubst $(HOLDERS_DIR)/%,%,$(dir $@)))..."
	@mkdir -p $(@D)
	@cp $(DATE_DIR)/$(notdir $(patsubst $(HOLDERS_DIR)/%,%,$(dir $@)))/outputs/ObserveOrder_$(DATE).txt $@
	@echo "✅ ObserveOrder file copied for band $(notdir $(patsubst $(HOLDERS_DIR)/%,%,$(dir $@)))"

# Simple copy target (doesn't depend on workflow)
copy_only:
	@echo "📋 Copying ObserveOrder files to holders directories..."
	@for band in $(BANDS); do \
		echo "📋 Copying ObserveOrder file for band $$band..."; \
		mkdir -p $(HOLDERS_DIR)/$$band/outputs; \
		cp $(DATE_DIR)/$$band/outputs/ObserveOrder_$(DATE).txt $(HOLDERS_DIR)/$$band/outputs/night_plan.csv; \
		echo "✅ ObserveOrder file copied for band $$band"; \
	done
	@echo "✅ All ObserveOrder files copied to holders directories!"

# Launch webapp
webapp:
	@echo "🌐 Launching AstroQ webapp..."
	@echo "📁 Using workdir as uptree path: $(WORKDIR)"
	@conda run -n astroq astroq kpfcc webapp -up $(WORKDIR)

# Clean up
clean:
	@echo " Cleaning up..."
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

.PHONY: all create_dirs clean status 