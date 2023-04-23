IMAGE_NAME=mobsuite
VERSION=3.1.4-mod

TAG1="us-central1-docker.pkg.dev/broad-dsp-lrma/general/$(IMAGE_NAME):$(VERSION)"
TAG2="us-central1-docker.pkg.dev/broad-dsp-lrma/general/$(IMAGE_NAME):latest"

all: build_and_push

build_and_push:
	docker buildx build --platform linux/amd64 --push -t $(TAG1) -t $(TAG2) .

build_and_push_nocache:
	docker buildx build --no-cache --platform linux/amd64 --push -t $(TAG1) -t $(TAG2) .
