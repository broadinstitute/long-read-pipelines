VERSION = 1.0
TAG1 = quay.io/ymostovoy/lr-differential-methylation:$(VERSION)
TAG2 = quay.io/ymostovoy/lr-differential-methylation:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
