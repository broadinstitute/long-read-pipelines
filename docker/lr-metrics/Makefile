TAG1 = quay.io/broad-long-read-pipelines/lr-metrics:0.01.07
TAG2 = quay.io/broad-long-read-pipelines/lr-metrics:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
