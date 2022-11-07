/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
/* angular unit corresponding to character */
    case 'h':
	msg("RA in hms, Dec in dms");		break;
    case 'd':
    case '°':
	msg("degrees");				break;
    case 'm':
    case '\'':
    case '´':
	msg("arcminutes");			break;
    case 's':
    case '"':
    case '¨':
	msg("arcseconds");			break;
    case 'r':
    default:
	msg("radians");				break;
