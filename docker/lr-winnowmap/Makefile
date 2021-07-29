VERSION = 2.03

TAG1 = us.gcr.io/broad-dsp-lrma/lr-winnowmap:$(VERSION)
# TAG2 = us.gcr.io/broad-dsp-lrma/lr-winnowmap:latest

all: build push

build:
	docker build \
	  -q \
	  -t $(TAG1) \
	  .

push:
	docker push $(TAG1)

