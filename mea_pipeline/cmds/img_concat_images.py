#!/usr/bin/env mea-pl

from mea_pipeline import cout, cerr, cexit, arg_parser


def init_argparser():
    p = arg_parser("Concatenate 2 or more images vertically/horizontally")

    p.add_argument("-o", "--outfile", default=None)
    p.add_argument("--layout", choices=["vertical", "horizontal"], default="vertical")
    p.add_argument("--resize", default=-1, type=float)
    p.add_argument("infiles", nargs="+")

    return p


def img_concat_images(args):

    from PIL import Image

    images = []
    for infile in args.infiles:
        cerr(f"Reading image: {infile}")
        images.append(Image.open(infile))

    new_img = None
    match args.layout:

        case "vertical":

            # get maximum width & total height
            max_width = max([img.width for img in images])
            total_height = sum([img.height for img in images])

            # generate new image
            new_img = Image.new("RGB", (max_width, total_height))
            off_y = 0
            for img in images:
                new_img.paste(img, (0, off_y))
                off_y += img.height

        case "horizontal":

            # get maximum height & total width
            max_height = max([img.height for img in images])
            total_width = sum([img.width for img in images])

            # generate new image
            new_img = Image.new("RGB", (total_width, max_height))
            off_x = 0
            for img in images:
                new_img.paste(img, (off_x, 0))
                off_x += img.width

        case _:
            raise ValueError("layout unknown")

    if args.resize > 0:
        cerr(f"Resize to {args.resize} scale.")
        new_img = new_img.resize(
            (int(new_img.width * args.resize), int(new_img.height * args.resize)),
            resample=Image.BICUBIC,
        )

    if new_img is not None and args.outfile:
        new_img.save(args.outfile)
        cerr(f"Writing image to: {args.outfile}")


def main(args):
    img_concat_images(args)


# EOF
