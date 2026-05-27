from case import build_case

if __name__ == "__main__":
    case, _ = build_case()
    case.clean()
    case.run()
