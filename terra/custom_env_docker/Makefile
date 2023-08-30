VERSION = 0.0.1
TAG1 = us.gcr.io/broad-dsp-lrma/lr-terra-base:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-terra-base:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
