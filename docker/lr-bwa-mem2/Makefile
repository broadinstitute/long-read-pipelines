VERSION = 2.2.1 # This should match the actual bwa-mem2 version
TAG1 = quay.io/broad-long-read-pipelines/lr-bwa-mem2:$(VERSION)
TAG2 = quay.io/broad-long-read-pipelines/lr-bwa-mem2:latest
TAG3 = us.gcr.io/broad-dsp-lrma/lr-bwa-mem2:$(VERSION)
TAG4 = us.gcr.io/broad-dsp-lrma/lr-bwa-mem2:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) -t $(TAG3) -t $(TAG4) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
	docker push $(TAG3)
	docker push $(TAG4)
