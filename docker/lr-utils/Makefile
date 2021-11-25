VERSION = 0.1.10
TAG1 = us.gcr.io/broad-dsp-lrma/lr-utils:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-utils:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
