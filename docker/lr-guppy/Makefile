IMAGE_NAME = lr-guppy
VERSION = 6.4.6 # This version number should follow the Guppy version number

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
