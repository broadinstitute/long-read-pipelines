TAG1 = quay.io/broad-long-read-pipelines/lr-asm-eval:0.01.02
TAG2 = quay.io/broad-long-read-pipelines/lr-asm-eval:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
