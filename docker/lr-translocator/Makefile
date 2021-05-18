VERSION = 0.1
#TAG1 = quay.io/broad-long-read-pipelines/lr-translocator:$(VERSION)
#TAG2 = quay.io/broad-long-read-pipelines/lr-translocator:latest
TAG3 = us.gcr.io/broad-dsp-lrma/lr-translocator:$(VERSION)
TAG4 = us.gcr.io/broad-dsp-lrma/lr-translocator:latest
#TAG1 = quay.io/ymostovoy/lr-translocator:$(VERSION)
#TAG2 = quay.io/ymostovoy/lr-translocator:latest


all: build push

build:
	docker build -t $(TAG3) -t $(TAG4) .

push:
	docker push $(TAG3)
	docker push $(TAG4)
