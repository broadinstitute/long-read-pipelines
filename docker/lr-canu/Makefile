
VERSION = 0.2.0

TAG1 = us.gcr.io/broad-dsp-lrma/lr-canu:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-canu:latest


all: | build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
