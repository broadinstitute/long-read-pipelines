VERSION = 2.2.1

TAG1 = us.gcr.io/broad-dsp-lrma/lr-sniffles2:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-sniffles2:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
