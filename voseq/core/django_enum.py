from abc import ABC, abstractmethod
from enum import unique, Enum
from types import MethodType
from typing import cast, Any, List, Tuple


class SupportsChoices(ABC):
    """ Class to add support for 'Interface like' classes having a get_choices method.
    Took as reference mypy's SupportsInt, for example.
    """
    @abstractmethod
    def get_choices(self) -> List[Tuple[Any, Any]]:
        ...


def django_enum(cls: Any) -> SupportsChoices:
    """ Allows us to use IntEnums as django choices
    IntEnum otherwise do not work natively.
    IntEnums cannot inherent, hence needs to be done with class decorator
    :param cls: the IntEnum class
    :return: the updated IntEnum class that can be used within django models
    """
    # ensure each of the enums are unique
    cls = unique(cls)
    # we need this to enable enums functioning in django templates
    cls.do_not_call_in_templates = True

    def __str__(self: Enum) -> str:
        # override this if you want to have other string representations than the variable name
        return str(self.value).replace('_', ' ')

    # we override the string representation to remove underscores from variable names
    cls.__str__ = __str__

    # add a convenience function to get django conform choices
    def get_choices(self: Any) -> List[Tuple[Any, Any]]:
        return [(choice.value, choice.name) for choice in self]

    cls.get_choices = MethodType(get_choices, cls)
    return cast(SupportsChoices, cls)
