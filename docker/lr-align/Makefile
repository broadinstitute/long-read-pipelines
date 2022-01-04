VERSION = 0.1.28
TAG1 = quay.io/broad-long-read-pipelines/lr-align:$(VERSION)
TAG2 = quay.io/broad-long-read-pipelines/lr-align:latest
TAG3 = us.gcr.io/broad-dsp-lrma/lr-align:$(VERSION)
TAG4 = us.gcr.io/broad-dsp-lrma/lr-align:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) -t $(TAG3) -t $(TAG4) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
	docker push $(TAG3)
	docker push $(TAG4)
