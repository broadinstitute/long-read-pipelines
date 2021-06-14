VERSION = 0.13
TAG1 = quay.io/broad-long-read-pipelines/lr-hifiasm:$(VERSION)
TAG2 = quay.io/broad-long-read-pipelines/lr-hifiasm:latest
TAG3 = us.gcr.io/broad-dsp-lrma/lr-hifiasm:$(VERSION)
TAG4 = us.gcr.io/broad-dsp-lrma/lr-hifiasm:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) -t $(TAG3) -t $(TAG4) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
	docker push $(TAG3)
	docker push $(TAG4) 
