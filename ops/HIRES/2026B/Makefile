# HIRES 2026B production Makefile. Run from repo root, e.g.:
#   export CC_OUTPUT_PATH=/path/to/output
#   conda activate astroq
#   make -C ops/HIRES/2026B DATE=2026-08-01 BANDS=full-band1

SHELL := /bin/bash

CONDA_ENV ?= astroq
WORKDIR ?= $(CC_OUTPUT_PATH)
SEMESTER_DIR = $(WORKDIR)/$(SEMESTER)
DATE_DIR = $(SEMESTER_DIR)/$(DATE)

MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
ENV_FILE ?= $(MAKEFILE_DIR).env
SHARED_ENV_FILE ?= $(MAKEFILE_DIR)../.env
LOAD_ENV = if [ -f "$(SHARED_ENV_FILE)" ]; then set -a; . "$(SHARED_ENV_FILE)"; set +a; fi; if [ -f "$(ENV_FILE)" ]; then set -a; . "$(ENV_FILE)"; set +a; fi

ALLOCATION_FILE ?= $(MAKEFILE_DIR)allocation_hires_cps.csv
REQUEST_URLS_CSV ?= $(MAKEFILE_DIR)request_urls.csv
CONFIG_TEMPLATE ?= $(MAKEFILE_DIR)../common/config_template.ini

SEMESTER ?= 2026B
START_DATE ?= 2026-08-01
END_DATE ?= 2027-01-31
BANDS ?= band1 full-band1
FILLER_PROGRAM ?= 2026B_E475

ifdef BAND
  BANDS := $(BAND)
endif

BOOST_UIDS ?=
BOOST_FACTOR ?= 5.0
ifneq ($(strip $(BOOST_UIDS)),)
  BOOST_ARGS = --boost "$(BOOST_UIDS)" $(BOOST_FACTOR)
else
  BOOST_ARGS =
endif

UNAME_S := $(shell uname)
SED_I_FLAG :=
ifeq ($(UNAME_S),Darwin)
  SED_I_FLAG := ''
endif

DATE ?= $(shell date +%Y-%m-%d)
ifdef date
	DATE := $(date)
endif

.PRECIOUS: $(foreach band,$(BANDS),$(DATE_DIR)/$(band)/config.ini)

all: $(foreach band,$(BANDS),$(DATE_DIR)/$(band)/plan-night-complete)
	@echo "Complete workflow finished for all bands!"
	@echo "Results available in:"
	@for band in $(BANDS); do \
		echo "  $(DATE_DIR)/$$band/"; \
	done
	@$(MAKE) combine_logs
	@$(MAKE) check_night_plans

$(DATE_DIR)/%/plan-night-complete: $(DATE_DIR)/%/plan-night-run
	@set -euo pipefail; \
	echo "Band $(notdir $(@D)) workflow complete — archiving webapp..."; \
	cd $(@D) && $(LOAD_ENV) && PYTHONUNBUFFERED=1 conda run --no-capture-output -n $(CONDA_ENV) astroq archive -cf config.ini 2>&1 | tee -a astroq.log; \
	touch $@

$(DATE_DIR)/%/plan-night-run: $(DATE_DIR)/%/plan-semester-run
	@set -euo pipefail; \
	echo "Running plan-night for band $(notdir $(@D))..."; \
	cd $(@D) && $(LOAD_ENV) && PYTHONUNBUFFERED=1 conda run --no-capture-output -n $(CONDA_ENV) astroq plan-night -cf config.ini 2>&1 | tee -a astroq.log; \
	touch $@

$(DATE_DIR)/%/plan-semester-run: $(DATE_DIR)/%/prep-run
	@set -euo pipefail; \
	echo "Running plan-semester for band $(notdir $(@D))..."; \
ifneq ($(strip $(BOOST_UIDS)),)
	echo "   Boost: $(BOOST_UIDS) (factor $(BOOST_FACTOR))"; \
endif
	cd $(@D) && $(LOAD_ENV) && PYTHONUNBUFFERED=1 conda run --no-capture-output -n $(CONDA_ENV) astroq plan-semester -cf config.ini $(BOOST_ARGS) 2>&1 | tee -a astroq.log; \
	touch $@

$(DATE_DIR)/%/prep-run: $(DATE_DIR)/%/config.ini
	@set -euo pipefail; \
	echo "Running prep for band $(notdir $(@D))..."; \
	test -f "$(REQUEST_URLS_CSV)"; \
	: > $(@D)/astroq.log; \
	QUEUE=$$(sed -n 's/^[[:space:]]*queue[[:space:]]*=[[:space:]]*//p' $(@D)/config.ini | tr -d ' \r\n'); \
	BAND_NUM=$$(echo $(notdir $(@D)) | sed 's/band//' | sed 's/full-//'); \
	IS_FULL_BAND=$$(echo $(notdir $(@D)) | grep -q '^full-' && echo "true" || echo "false"); \
	if [ "$$IS_FULL_BAND" = "true" ]; then \
		echo "Prep full-band $$BAND_NUM (queue=$$QUEUE)..."; \
		cd $(@D) && $(LOAD_ENV) && PYTHONUNBUFFERED=1 conda run --no-capture-output -n $(CONDA_ENV) astroq prep $$QUEUE -cf config.ini -ru "$(REQUEST_URLS_CSV)" -fillers $(FILLER_PROGRAM) -band $$BAND_NUM -full 2>&1 | tee -a astroq.log; \
	else \
		echo "Prep band $$BAND_NUM (queue=$$QUEUE)..."; \
		cd $(@D) && $(LOAD_ENV) && PYTHONUNBUFFERED=1 conda run --no-capture-output -n $(CONDA_ENV) astroq prep $$QUEUE -cf config.ini -ru "$(REQUEST_URLS_CSV)" -fillers $(FILLER_PROGRAM) -band $$BAND_NUM 2>&1 | tee -a astroq.log; \
	fi; \
	touch $@

$(DATE_DIR)/%/config.ini: create_dirs
	@echo "Copying config template for band $(notdir $(@D))..."
	@mkdir -p $(@D)
	@cp $(CONFIG_TEMPLATE) $@
	@sed -i $(SED_I_FLAG) "s|CURRENT_DATE_PLACEHOLDER|$(DATE)|g" $@
	@sed -i $(SED_I_FLAG) "s|START_DATE_PLACEHOLDER|$(START_DATE)|g" $@
	@sed -i $(SED_I_FLAG) "s|END_DATE_PLACEHOLDER|$(END_DATE)|g" $@
	@sed -i $(SED_I_FLAG) "s|SEMESTER_PLACEHOLDER|$(SEMESTER)|g" $@
	@sed -i $(SED_I_FLAG) "s|WORKDIR_PLACEHOLDER|$(@D)|g" $@
	@echo "Config file created for band $(notdir $(@D))"

validate_date:
ifdef date
	@if ! echo "$(date)" | grep -E '^[0-9]{4}-[0-9]{2}-[0-9]{2}$$' > /dev/null; then \
		echo "Error: Invalid date format: $(date) (use YYYY-MM-DD)"; exit 1; \
	fi
	@if ! python3 -c "from datetime import datetime; datetime.strptime('$(date)', '%Y-%m-%d')" > /dev/null 2>&1; then \
		echo "Error: Invalid date: $(date)"; exit 1; \
	fi
else
	@echo "Using date: $(DATE)"
endif

create_dirs: validate_date
	@if [ -z "$(strip $(WORKDIR))" ]; then \
		echo "Error: set CC_OUTPUT_PATH or WORKDIR to the output parent directory"; exit 1; \
	fi
	@mkdir -p $(DATE_DIR)

combine_logs:
	@echo "# Combined AstroQ Logs for $(DATE)" > $(DATE_DIR)/all_band_logs.log
	@echo "# Bands: $(BANDS)" >> $(DATE_DIR)/all_band_logs.log
	@echo "" >> $(DATE_DIR)/all_band_logs.log
	@for band in $(BANDS); do \
		if [ -f "$(DATE_DIR)/$$band/astroq.log" ]; then \
			echo "## Band: $$band" >> $(DATE_DIR)/all_band_logs.log; \
			cat $(DATE_DIR)/$$band/astroq.log >> $(DATE_DIR)/all_band_logs.log; \
			echo "" >> $(DATE_DIR)/all_band_logs.log; \
		fi; \
	done

check_night_plans:
	@set -euo pipefail; \
	$(LOAD_ENV) && PYTHONUNBUFFERED=1 conda run --no-capture-output -n $(CONDA_ENV) python -m astroq.scripts.check_night_plans -s $(SEMESTER) -d $(DATE) $(DATE_DIR)

webapp:
	@chmod +x $(WORKDIR)/launch_webapp.sh
	@$(WORKDIR)/launch_webapp.sh

clean:
	@echo "This will remove: $(DATE_DIR)"
	@echo -n "Are you sure? (y/N): " && read -r confirm && [ "$$confirm" = "y" ] || exit 1
	@rm -rf $(DATE_DIR)

status:
	@echo "Semester:     $(SEMESTER)"
	@echo "Date:         $(DATE)"
	@echo "Start/End:    $(START_DATE) .. $(END_DATE)"
	@echo "Bands:        $(BANDS)"
	@echo "Filler:       $(FILLER_PROGRAM)"
	@echo "Conda env:    $(CONDA_ENV)"
	@echo "Request URLs: $(REQUEST_URLS_CSV)"
	@echo "WORKDIR:      $(WORKDIR)"
	@echo "DATE_DIR:     $(DATE_DIR)"

complete: all webapp

.PHONY: all create_dirs clean status webapp complete combine_logs check_night_plans prep plan-semester plan-night

prep: $(foreach band,$(BANDS),$(DATE_DIR)/$(band)/prep-run)
plan-semester: $(foreach band,$(BANDS),$(DATE_DIR)/$(band)/plan-semester-run)
plan-night: $(foreach band,$(BANDS),$(DATE_DIR)/$(band)/plan-night-run)
