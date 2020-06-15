VERSION = 0.1.0
TAG1 = quay.io/broad-long-read-pipelines/lr-softclipper:$(VERSION)
TAG2 = quay.io/broad-long-read-pipelines/lr-softclipper:latest
TAG3 = us.gcr.io/broad-dsp-lrma/lr-softclipper:$(VERSION)
TAG4 = us.gcr.io/broad-dsp-lrma/lr-softclipper:latest

all: build push

build:
	cp ../../scripts/soft_clipper.py .
	docker build -t $(TAG1) -t $(TAG2) -t $(TAG3) -t $(TAG4) .
	rm soft_clipper.py

push:
	docker push $(TAG1)
	docker push $(TAG2)
	docker push $(TAG3)
	docker push $(TAG4)
