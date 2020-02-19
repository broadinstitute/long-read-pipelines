TAG1 = quay.io/broad-long-read-pipelines/sr-utils:0.01.00
TAG2 = quay.io/broad-long-read-pipelines/sr-utils:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
