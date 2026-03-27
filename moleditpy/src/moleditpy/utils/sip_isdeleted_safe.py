try:
    from PyQt6 import sip as _sip
    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except ImportError:
    try:
        import sip as _sip
        _sip_isdeleted = getattr(_sip, "isdeleted", None)
    except ImportError:
        _sip = None
        _sip_isdeleted = None


def sip_isdeleted_safe(obj) -> bool:
    """Return True if sip reports the given wrapper object as deleted.

    This function is conservative: if SIP isn't available or any error
    occurs while checking, it returns False (i.e. not deleted) so that the
    caller can continue other lightweight guards (like checking scene()).
    """
    try:
        if obj is None:
            return True
        if _sip_isdeleted is None:
            return False
        return bool(_sip_isdeleted(obj))
    except (AttributeError, RuntimeError, TypeError):
        return False
