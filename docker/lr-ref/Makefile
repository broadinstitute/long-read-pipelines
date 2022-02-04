VERSION = 0.1.0
TAG1 = us.gcr.io/broad-dsp-lrma/lr-ref:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-ref:latest

all: | build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

build_no_cache:
	docker build --no-cache -t $(TAG1) -t $(TAG2) .


push:
	docker push $(TAG1)
	docker push $(TAG2)
