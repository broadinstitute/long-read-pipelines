VERSION = 2.26-gcloud
TAG1 = us.gcr.io/broad-dsp-lrma/lr-minimap2:$(VERSION)

all: build push

build:
	docker build -t $(TAG1) .

push:
	docker push $(TAG1)
