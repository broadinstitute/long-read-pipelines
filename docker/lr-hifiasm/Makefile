IMAGE_NAME = lr-hifiasm
VERSION = 0.20.0

TAG1 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):latest

all: | build push

build:
	docker buildx build --platform linux/amd64            -t $(TAG1) -t $(TAG2) .

build_no_cache:
	docker buildx build --platform linux/amd64 --no-cache -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
