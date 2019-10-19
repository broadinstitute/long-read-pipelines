TAG1 = quay.io/broad-long-read-pipelines/lr-asm:0.01.10
TAG2 = quay.io/broad-long-read-pipelines/lr-asm:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
