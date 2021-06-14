VERSION = 0.3.1 # This should match the MosDepth version number
TAG1 = quay.io/broad-long-read-pipelines/lr-mosdepth:$(VERSION)
TAG2 = quay.io/broad-long-read-pipelines/lr-mosdepth:latest
TAG3 = us.gcr.io/broad-dsp-lrma/lr-mosdepth:$(VERSION)
TAG4 = us.gcr.io/broad-dsp-lrma/lr-mosdepth:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) -t $(TAG3) -t $(TAG4) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
	docker push $(TAG3)
	docker push $(TAG4)
