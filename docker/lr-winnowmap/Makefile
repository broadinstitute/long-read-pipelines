VERSION = 2.03
#TAG1 = quay.io/broad-long-read-pipelines/lr-winnowmap:$(VERSION)
#TAG2 = quay.io/broad-long-read-pipelines/lr-winnowmap:latest
TAG3 = us.gcr.io/broad-dsp-lrma/lr-winnowmap:$(VERSION)
TAG4 = us.gcr.io/broad-dsp-lrma/lr-winnowmap:latest
#TAG1 = quay.io/ymostovoy/lr-winnowmap:$(VERSION)
#TAG2 = quay.io/ymostovoy/lr-winnowmap:latest

all: build push

build:
	docker build -t $(TAG3) -t $(TAG4) .

push:
	docker push $(TAG3)
	docker push $(TAG4)
