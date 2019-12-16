
# Getting Cosmic Data programmatically
if args.cosmic_user and args.cosmic_pass:
    r = requests.get(COSMIC_URL, auth=(args.cosmic_user, args.cosmic_pass))
    cosmic_data = requests.get(r.json()['url'], stream=True)