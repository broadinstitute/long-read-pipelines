VERSION = 0.1.3
TAG1 = us.gcr.io/broad-dsp-lrma/lr-basic:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-basic:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
