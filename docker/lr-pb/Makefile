TAG1 = quay.io/broad-long-read-pipelines/lr-pb:0.01.04
TAG2 = quay.io/broad-long-read-pipelines/lr-pb:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
