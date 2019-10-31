TAG1 = quay.io/broad-long-read-pipelines/lr-sv:0.1.0
TAG2 = quay.io/broad-long-read-pipelines/lr-sv:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
