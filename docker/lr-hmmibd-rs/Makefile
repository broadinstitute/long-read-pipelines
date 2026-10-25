IMAGE_NAME = lr-hmmibd-rs
# Version is pinned to the hmmibd-rs crate version: the image is built from hmmibd-rs tag v$(VERSION)
# (passed as --build-arg below) and tagged :$(VERSION). It also carries the bcf2hmmibd converter and
# bcftools/tabix/gawk so filtration + conversion run in-image.
VERSION = 0.1.7

TAG1 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):latest

all: | build push

build:
	docker build --build-arg HMMIBD_RS_REF=v$(VERSION) -t $(TAG1) -t $(TAG2) .

build_no_cache:
	docker build --no-cache --build-arg HMMIBD_RS_REF=v$(VERSION) -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
