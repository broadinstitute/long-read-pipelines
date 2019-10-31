TAG1 = quay.io/broad-long-read-pipelines/lr-ont:0.01.01
TAG2 = quay.io/broad-long-read-pipelines/lr-ont:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
