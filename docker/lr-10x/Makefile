TAG1 = quay.io/broad-long-read-pipelines/lr-10x:0.01.08
TAG2 = quay.io/broad-long-read-pipelines/lr-10x:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
