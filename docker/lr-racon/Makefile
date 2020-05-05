IMAGE_NAME = lr-racon
VERSION = 0.1.0

TAG1 = quay.io/broad-long-read-pipelines/$(IMAGE_NAME):$(VERSION)
TAG2 = quay.io/broad-long-read-pipelines/$(IMAGE_NAME):latest
TAG3 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):$(VERSION)
TAG4 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):latest

all: | build push

build:
	docker build -t $(TAG1) -t $(TAG2) -t $(TAG3) -t $(TAG4) .

build_no_cache:
	docker build --no-cache -t $(TAG1) -t $(TAG2) -t $(TAG3) -t $(TAG4) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
	docker push $(TAG3)
	docker push $(TAG4)
