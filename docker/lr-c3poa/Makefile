TAG1 = quay.io/broad-long-read-pipelines/lr-c3poa:0.01.03
TAG2 = quay.io/broad-long-read-pipelines/lr-c3poa:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
