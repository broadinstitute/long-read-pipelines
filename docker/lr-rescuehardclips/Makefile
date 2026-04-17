VERSION = 0.1.0 
TAG1 = us.gcr.io/broad-dsp-lrma/lr-rescuehardclips:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-rescuehardclips:latest

ARCH := $(shell uname -m)
BUILDER_NAME := lr-builder

# Set platform flag only when needed
ifeq ($(ARCH), arm64)
	PLATFORM_FLAG := --platform linux/amd64
else
	PLATFORM_FLAG :=
endif

.PHONY: all build push

all: build push

build:
	@docker buildx inspect $(BUILDER_NAME) > /dev/null 2>&1 || \
		docker buildx create --name $(BUILDER_NAME) --use
	@docker buildx use $(BUILDER_NAME)
	docker buildx build $(PLATFORM_FLAG) \
		--load \
		-t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
