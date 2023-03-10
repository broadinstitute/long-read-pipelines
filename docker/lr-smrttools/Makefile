VERSION = 12.0.0.176214
TAG1 = us.gcr.io/broad-dsp-lrma/lr-smrttools:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-smrttools:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
