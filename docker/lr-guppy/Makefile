VERSION = 6.0.1 # This version number should follow the Guppy version number
TAG1 = us.gcr.io/broad-dsp-lrma/lr-guppy:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-guppy:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
