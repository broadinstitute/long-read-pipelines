TAG1 = quay.io/broad-long-read-pipelines/lr-utils:0.01.04
TAG2 = quay.io/broad-long-read-pipelines/lr-utils:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
