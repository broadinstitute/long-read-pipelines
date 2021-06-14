VERSION = 4.5.2 # This version number should follow the Guppy version number
TAG1 = quay.io/broad-long-read-pipelines/lr-guppy:$(VERSION)
TAG2 = quay.io/broad-long-read-pipelines/lr-guppy:latest
TAG3 = us.gcr.io/broad-dsp-lrma/lr-guppy:$(VERSION)
TAG4 = us.gcr.io/broad-dsp-lrma/lr-guppy:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) -t $(TAG3) -t $(TAG4) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
	docker push $(TAG3)
	docker push $(TAG4)
