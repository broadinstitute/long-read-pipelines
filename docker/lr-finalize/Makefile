TAG1 = quay.io/broad-long-read-pipelines/lr-finalize:0.01.00
TAG2 = quay.io/broad-long-read-pipelines/lr-finalize:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
