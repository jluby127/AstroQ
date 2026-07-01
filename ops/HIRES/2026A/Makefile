# HIRES 2026A production pipeline. Run from repo root:
#   make -C ops/HIRES/2026A DATE=2026-07-03 BANDS=band1

ASTROQ_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../../..)
SEMESTER ?= 2026A
DATE ?=
BANDS ?= band1
UPTREE ?= $(CC_OUTPUT_PATH)
FILLER ?= 2026A_E473
REQUEST_URLS ?= ops/HIRES/2026A/request_urls.csv
START_DATE ?= 2026-02-01
END_DATE ?= 2026-07-31

ifeq ($(strip $(UPTREE)),)
$(error Set CC_OUTPUT_PATH or UPTREE to the output parent directory)
endif

ifeq ($(strip $(DATE)),)
$(error DATE is required, e.g. make DATE=2026-07-03)
endif

BAND := $(firstword $(BANDS))
BAND_NUM := $(patsubst band%,%,$(BAND))
RUNDIR := $(abspath $(UPTREE)/$(SEMESTER)/$(DATE)/$(BAND))
CONFIG := $(RUNDIR)/config.ini
TEMPLATE := $(ASTROQ_ROOT)/ops/HIRES/common/config_template.ini

.PHONY: all prep plan-semester plan-night archive config

all: archive

config: $(CONFIG)

$(CONFIG):
	mkdir -p $(RUNDIR)
	sed -e 's|WORKDIR_PLACEHOLDER|$(RUNDIR)|' \
	    -e 's|CURRENT_DATE_PLACEHOLDER|$(DATE)|' \
	    -e 's|START_DATE_PLACEHOLDER|$(START_DATE)|' \
	    -e 's|END_DATE_PLACEHOLDER|$(END_DATE)|' \
	    -e 's|SEMESTER_PLACEHOLDER|$(SEMESTER)|' \
	    $(TEMPLATE) > $(CONFIG)

prep: $(CONFIG)
	cd $(ASTROQ_ROOT) && astroq prep hirescps -cf $(CONFIG) -band $(BAND_NUM) \
	    -ru $(REQUEST_URLS) -fillers $(FILLER)

plan-semester: prep
	cd $(ASTROQ_ROOT) && astroq plan-semester -cf $(CONFIG)

plan-night: plan-semester
	cd $(ASTROQ_ROOT) && astroq plan-night -cf $(CONFIG)

archive: plan-night
	cd $(ASTROQ_ROOT) && astroq archive -cf $(CONFIG)
