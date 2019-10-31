TAG1 = quay.io/broad-long-read-pipelines/lr-align:0.01.24
TAG2 = quay.io/broad-long-read-pipelines/lr-align:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
