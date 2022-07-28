IMAGE_NAME = lr-jupyter
VERSION = 0.0.3

INTERACTIVE_SUFFIX = _interactive

TAG1 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):latest

TAG3 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME)$(INTERACTIVE_SUFFIX):$(VERSION)
TAG4 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME)$(INTERACTIVE_SUFFIX):latest

all: | build push

build:
	docker build -t $(TAG1) -t $(TAG2)  .
	cat Dockerfile | sed 's#CMD \["start-notebook.sh"\]#CMD \["/bin/bash"\]#g' > Dockerfile.interactive
	docker build - -t $(TAG3) -t $(TAG4) . < Dockerfile.interactive

build_no_cache:
	docker build --no-cache -t $(TAG1) -t $(TAG2) .
	cat Dockerfile | sed 's#CMD \["start-notebook.sh"\]#CMD \["/bin/bash"\]#g' > Dockerfile.interactive
	docker build - -t $(TAG3) -t $(TAG4) . < Dockerfile.interactive

push:
	docker push $(TAG1)
	docker push $(TAG2)
	docker push $(TAG3)
	docker push $(TAG4)
