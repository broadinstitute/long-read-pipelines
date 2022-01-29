VERSION = 1.2.1 # This should match the WhatsHap version number
TAG1 = us.gcr.io/broad-dsp-lrma/lr-whatshap:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-whatshap:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
