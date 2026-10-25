IMAGE_NAME = lr-hmmibd-rs
# Image revision. Tracks the hmmibd-rs crate version (0.1.7) with a trailing image rev: this build
# also carries the bcf2hmmibd converter and bcftools/tabix/gawk so filtration+conversion run in-image.
VERSION = 0.1.8

TAG1 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):latest

all: | build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

build_no_cache:
	docker build --no-cache -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
