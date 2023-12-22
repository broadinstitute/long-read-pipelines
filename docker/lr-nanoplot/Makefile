IMAGE_NAME = lr-nanoplot
VERSION = 1.42.0
TAG1 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
