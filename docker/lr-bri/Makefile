TAG1 = quay.io/broad-long-read-pipelines/lr-bri:0.1.17
TAG2 = quay.io/broad-long-read-pipelines/lr-bri:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
