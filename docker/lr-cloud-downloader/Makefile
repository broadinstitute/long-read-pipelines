TAG1 = kgarimella/lr-cloud-downloader:0.02.0
TAG2 = kgarimella/lr-cloud-downloader:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
