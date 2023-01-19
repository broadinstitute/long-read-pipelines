VERSION = 1.3.0 # This should match the version number at https://dalexander.github.io/admixture/download.html

TAG1 = us.gcr.io/broad-dsp-lrma/lr-admixture:$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/lr-admixture:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
