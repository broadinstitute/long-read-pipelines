VERSION = 0.1.34
TAG1 = us.gcr.io/broad-dsp-lrma/lr-pb:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-pb:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
