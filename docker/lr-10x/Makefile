VERSION = 0.1.16
NAME = lr-10x

TAG1 = quay.io/broad-long-read-pipelines/$(NAME):$(VERSION)
TAG2 = quay.io/broad-long-read-pipelines/$(NAME):latest
TAG3 = us.gcr.io/broad-dsp-lrma/$(NAME):$(VERSION)
TAG4 = us.gcr.io/broad-dsp-lrma/$(NAME):latest

all: | build push

build:
	docker build -t $(TAG1) -t $(TAG2) -t $(TAG3) -t $(TAG4) .

build_no_cache:
	docker build --no-cache -t $(TAG1) -t $(TAG2) -t $(TAG3) -t $(TAG4) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
	docker push $(TAG3)
	docker push $(TAG4)
