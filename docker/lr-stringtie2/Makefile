VERSION = 2.1.6 # This should match the Stringtie version number
TAG1 = quay.io/broad-long-read-pipelines/lr-stringtie2:$(VERSION)
TAG2 = quay.io/broad-long-read-pipelines/lr-stringtie2:latest
TAG3 = us.gcr.io/broad-dsp-lrma/lr-stringtie2:$(VERSION)
TAG4 = us.gcr.io/broad-dsp-lrma/lr-stringtie2:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) -t $(TAG3) -t $(TAG4) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
	docker push $(TAG3)
	docker push $(TAG4)
