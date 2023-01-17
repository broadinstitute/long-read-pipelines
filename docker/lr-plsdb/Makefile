IMAGE_NAME=plsdb
VERSION=0.1.1

TAG1="us-east1-docker.pkg.dev/broad-dsp-lrma/fusilli/$(IMAGE_NAME):$(VERSION)"
TAG2="us-east1-docker.pkg.dev/broad-dsp-lrma/fusilli/$(IMAGE_NAME):latest"

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

build_no_cache:
	docker build --no-cache -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
